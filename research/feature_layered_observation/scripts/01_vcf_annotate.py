#!/usr/bin/env python3
"""
Phase B1 — VCF annotation for merged master.

Purpose:
  Attach ClairS / ClairS-TO caller signals (INFO + FORMAT tags) to every row of
  the Thread B merged master (7 samples x {paired_full, to_pileup}).

Input:
  - Merged master tsv.gz (ISM aggregate features only; key = sample, mode, Chr, Pos)

Output:
  - merged_with_vcf.tsv.gz under research/feature_layered_observation/data/
    (original columns + vcf_* columns; row count preserved)
  - per-sample join stats (printed + logged)
  - hcc1395 caller_AF vs ISM AF Pearson correlation on TP rows

Strategy:
  - Per (sample, mode), iterate both filtered_snv_tp.vcf.gz and filtered_snv_fp.vcf.gz.
  - For each variant, build a record key = (sample, mode, chr, pos).
  - Harmonize INFO and FORMAT tags across ClairS (paired) vs ClairS-TO schemas.
  - Use multiprocessing.Pool(14) across (sample, mode) pairs.
  - After all VCF frames collected, merge into master via left join.

Design notes:
  - paired_full VCF uses ClairS v0.4.0 schema. Key INFO: H flag, FAU..RTU counts.
    Key FORMAT: GT, GQ, DP, AF, AD, NAF, NDP, NAD, AU, CU, GU, TU, NAU..NTU.
  - to_pileup VCF uses ClairS-TO v0.3.0 schema. Key INFO: H, FAU..RTU, SB,
    PoN_1..PoN_4, Verdict_Germline, Verdict_Somatic, Verdict_SubclonalSomatic.
    Key FORMAT: GT, GQ, DP, AF, AD, AU..TU (no NAF/NDP/NAD since no normal BAM).
  - Strand bias (SB): only present in ClairS-TO (paired has none in INFO).
  - If a key appears in both TP and FP VCFs (should not, but guard), last wins
    and a warning is logged.
"""

from __future__ import annotations

import gzip
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pysam

ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
MASTER = ROOT / "research/ng_kde_rescaling/data/merged_7samples_paired_full_plus_7to_full_20260423.tsv.gz"
OUT_DIR = ROOT / "research/feature_layered_observation/data"
LOG_DIR = ROOT / "research/feature_layered_observation/logs"
OUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# VCF path table per (sample, mode). Keys: filtered_snv_tp.vcf.gz, filtered_snv_fp.vcf.gz.
# mode values follow merged master: paired_full, to_pileup.
# ---------------------------------------------------------------------------
def _paired_full(dir_name: str) -> Path:
    # canonical/<sample>/paired_full/<dir>/longphase_s/filtered_snv_{tp,fp}.vcf.gz
    return ROOT / "output/canonical" / dir_name

VCF_PATHS: dict[tuple[str, str], dict[str, Path]] = {
    # --- paired_full (ClairS v0.4.0) ---
    ("HCC1395", "paired_full"): {
        "dir": ROOT / "output/canonical/HCC1395/paired_full/20260420_HCC1395_paired_full_full_2/longphase_s",
    },
    ("HCC1395_DORADO", "paired_full"): {
        "dir": ROOT / "output/canonical/HCC1395_DORADO/paired_full/20260420_HCC1395_DORADO_paired_full_full/longphase_s",
    },
    ("COLO829", "paired_full"): {
        "dir": ROOT / "output/canonical/COLO829/paired_full/20260421_COLO829_paired_full_full/longphase_s",
    },
    ("H2009", "paired_full"): {
        "dir": ROOT / "output/canonical/H2009/paired_full/20260421_H2009_paired_full_full/longphase_s",
    },
    ("H1437", "paired_full"): {
        "dir": ROOT / "output/canonical/H1437/paired_full/20260421_H1437_paired_full_full/longphase_s",
    },
    ("HCC1937", "paired_full"): {
        "dir": ROOT / "output/canonical/HCC1937/paired_full/20260421_HCC1937_paired_full_full/longphase_s",
    },
    ("HCC1954", "paired_full"): {
        "dir": ROOT / "output/canonical/HCC1954/paired_full/20260421_HCC1954_paired_full_full/longphase_s",
    },

    # --- to_pileup (ClairS-TO v0.3.0) ---
    ("HCC1395", "to_pileup"): {
        "dir": ROOT / "output/synthesis/research_rounds/20260315_hcc1395_to_pilot/step04_benchmark_longphase_to",
    },
    ("HCC1395_DORADO", "to_pileup"): {
        "dir": ROOT / "output/synthesis/research_rounds/archive/202603_early_pilots/20260315_hcc1395_dorado_to_pilot/step04_benchmark_longphase_to",
    },
    ("COLO829", "to_pileup"): {
        "dir": ROOT / "output/synthesis/research_rounds/20260423_colo829_to_pilot/step04_benchmark_longphase_to",
    },
    ("H2009", "to_pileup"): {
        "dir": ROOT / "output/synthesis/research_rounds/archive/202603_early_pilots/20260318_h2009_to_pilot_fastresume/step04_benchmark_longphase_to",
    },
    ("H1437", "to_pileup"): {
        "dir": ROOT / "output/synthesis/research_rounds/archive/202603_early_pilots/20260318_h1437_to_pilot_fastresume/step04_benchmark_longphase_to",
    },
    ("HCC1937", "to_pileup"): {
        "dir": ROOT / "output/synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1937_to_pilot_fastresume/step04_benchmark_longphase_to",
    },
    ("HCC1954", "to_pileup"): {
        "dir": ROOT / "output/synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1954_to_pilot_fastresume/step04_benchmark_longphase_to",
    },
}

# ---------------------------------------------------------------------------
# Tags extracted from VCF. Names normalized to "vcf_<TAG>" in final dataframe.
# Flags (0/1) produce numeric 0 or 1; missing tags become NaN.
# ---------------------------------------------------------------------------
INFO_FLAGS_TO = [
    "H",
    "PoN_1", "PoN_2", "PoN_3", "PoN_4",
    "Verdict_Germline", "Verdict_Somatic", "Verdict_SubclonalSomatic",
]
INFO_FLAGS_PAIRED = ["H"]

INFO_NUMERIC_COMMON = [
    "FAU", "FCU", "FGU", "FTU",
    "RAU", "RCU", "RGU", "RTU",
]
INFO_NUMERIC_TO_ONLY = ["SB"]

FORMAT_COMMON = ["GT", "GQ", "DP", "AF"]
FORMAT_PAIRED_EXTRA = ["NAF", "NDP"]
# AD / NAD are tuples -> store as alt allelic depth + ref allelic depth separately.


def _flag_present(rec: pysam.VariantRecord, key: str) -> int | float:
    # Flags in pysam appear as key in rec.info with True value; absence means NaN/0 depending on design.
    try:
        val = rec.info.get(key, None)
    except KeyError:
        return 0
    if val is None:
        return 0
    # pysam returns True for flags, or a scalar; treat any truthy as 1.
    return 1 if bool(val) else 0


def _info_numeric(rec: pysam.VariantRecord, key: str) -> float:
    try:
        v = rec.info.get(key, None)
    except KeyError:
        return np.nan
    if v is None:
        return np.nan
    if isinstance(v, tuple):
        return float(v[0]) if len(v) > 0 and v[0] is not None else np.nan
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


def _format_scalar(sample_field: dict[str, Any], key: str) -> Any:
    v = sample_field.get(key, None)
    if v is None:
        return np.nan
    if isinstance(v, tuple):
        # e.g. AD/NAD -> tuple (ref, alt). Return alt count.
        return v
    return v


def _parse_record(rec: pysam.VariantRecord, mode: str, label: int) -> dict[str, Any]:
    """Turn one VCF record into a flat dict keyed by vcf_* column names."""
    out: dict[str, Any] = {}
    out["Chr"] = rec.chrom
    out["Pos"] = rec.pos  # 1-based, matches master Pos
    out["vcf_Ref"] = rec.ref
    out["vcf_Alt"] = rec.alts[0] if rec.alts else None
    out["vcf_QUAL"] = rec.qual if rec.qual is not None else np.nan
    out["vcf_FILTER"] = ";".join(rec.filter.keys()) if rec.filter.keys() else "PASS"
    out["vcf_tp_label_from_vcf"] = label  # 1 for TP, 0 for FP

    # Flags
    flags = INFO_FLAGS_TO if mode == "to_pileup" else INFO_FLAGS_PAIRED
    for flag in flags:
        out[f"vcf_{flag}"] = _flag_present(rec, flag)

    # Numeric INFO
    numeric_keys = INFO_NUMERIC_COMMON + (INFO_NUMERIC_TO_ONLY if mode == "to_pileup" else [])
    for k in numeric_keys:
        out[f"vcf_{k}"] = _info_numeric(rec, k)

    # FORMAT (single sample VCFs)
    if len(rec.samples) > 0:
        sname = list(rec.samples.keys())[0]
        sf = rec.samples[sname]
        for k in FORMAT_COMMON:
            val = _format_scalar(sf, k)
            if k == "GT":
                # Encode GT tuple as string
                if isinstance(val, tuple):
                    val = "|".join(["." if v is None else str(v) for v in val])
            out[f"vcf_{k}"] = val
        if mode == "paired_full":
            for k in FORMAT_PAIRED_EXTRA:
                out[f"vcf_{k}"] = _format_scalar(sf, k)
        # AD -> split into alt depth and ref depth
        ad = _format_scalar(sf, "AD")
        if isinstance(ad, tuple) and len(ad) >= 2:
            out["vcf_AD_ref"] = ad[0]
            out["vcf_AD_alt"] = ad[1]
        else:
            out["vcf_AD_ref"] = np.nan
            out["vcf_AD_alt"] = np.nan
        if mode == "paired_full":
            nad = _format_scalar(sf, "NAD")
            if isinstance(nad, tuple) and len(nad) >= 2:
                out["vcf_NAD_ref"] = nad[0]
                out["vcf_NAD_alt"] = nad[1]
            else:
                out["vcf_NAD_ref"] = np.nan
                out["vcf_NAD_alt"] = np.nan

    return out


def _parse_vcf(path: Path, mode: str, label: int) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    rows: list[dict[str, Any]] = []
    with pysam.VariantFile(str(path)) as vf:
        for rec in vf:
            rows.append(_parse_record(rec, mode, label))
    return rows


def _worker(args: tuple[str, str, dict[str, Path]]) -> tuple[str, str, pd.DataFrame, dict[str, Any]]:
    sample, mode, paths = args
    t0 = time.time()
    tp_path = paths["dir"] / "filtered_snv_tp.vcf.gz"
    fp_path = paths["dir"] / "filtered_snv_fp.vcf.gz"
    stats: dict[str, Any] = {
        "sample": sample,
        "mode": mode,
        "tp_vcf": str(tp_path),
        "fp_vcf": str(fp_path),
        "tp_exists": tp_path.exists(),
        "fp_exists": fp_path.exists(),
        "tp_rows": 0,
        "fp_rows": 0,
        "duplicate_keys": 0,
    }
    tp_rows = _parse_vcf(tp_path, mode, 1)
    fp_rows = _parse_vcf(fp_path, mode, 0)
    stats["tp_rows"] = len(tp_rows)
    stats["fp_rows"] = len(fp_rows)
    stats["elapsed_sec"] = round(time.time() - t0, 2)

    combined = tp_rows + fp_rows
    if not combined:
        return sample, mode, pd.DataFrame(), stats

    df = pd.DataFrame(combined)
    df["sample"] = sample
    df["mode"] = mode

    # Guard against duplicate (sample, mode, Chr, Pos) – prefer TP over FP.
    dup_mask = df.duplicated(subset=["sample", "mode", "Chr", "Pos"], keep=False)
    stats["duplicate_keys"] = int(dup_mask.sum())
    if stats["duplicate_keys"] > 0:
        # Sort so TP (label=1) comes first, then drop duplicates keeping first.
        df = df.sort_values("vcf_tp_label_from_vcf", ascending=False).drop_duplicates(
            subset=["sample", "mode", "Chr", "Pos"], keep="first"
        )
    return sample, mode, df, stats


def main() -> int:
    print(f"[{time.strftime('%H:%M:%S')}] Loading master: {MASTER}")
    master = pd.read_csv(MASTER, sep="\t")
    print(f"  master rows = {len(master):,}")
    print(f"  master cols = {list(master.columns)}")

    tasks = [(s, m, p) for (s, m), p in VCF_PATHS.items()]
    print(f"[{time.strftime('%H:%M:%S')}] Spawning {min(len(tasks), 14)} workers for {len(tasks)} (sample,mode) tasks")
    with mp.Pool(processes=min(len(tasks), 14)) as pool:
        results = pool.map(_worker, tasks)

    vcf_dfs: list[pd.DataFrame] = []
    stats_rows: list[dict[str, Any]] = []
    for sample, mode, df, stats in results:
        stats_rows.append(stats)
        flag = "OK" if not df.empty else "MISS"
        print(
            f"  [{flag}] {sample}/{mode}: TP={stats['tp_rows']}, FP={stats['fp_rows']}, "
            f"dup={stats['duplicate_keys']}, t={stats['elapsed_sec']}s"
        )
        if not df.empty:
            vcf_dfs.append(df)

    stats_df = pd.DataFrame(stats_rows)
    stats_path = LOG_DIR / "01_vcf_parse_stats.tsv"
    stats_df.to_csv(stats_path, sep="\t", index=False)
    print(f"[{time.strftime('%H:%M:%S')}] Wrote parse stats -> {stats_path}")

    if not vcf_dfs:
        print("ERROR: no VCF dataframes parsed – aborting.")
        return 1

    vcf_all = pd.concat(vcf_dfs, ignore_index=True, sort=False)
    print(f"  Concatenated VCF rows = {len(vcf_all):,}")

    # ---- Join to master ----
    # Master has columns (Chr, Pos) as plain strings/ints. Ensure dtypes align.
    master["Pos"] = master["Pos"].astype(int)
    vcf_all["Pos"] = vcf_all["Pos"].astype(int)
    master["Chr"] = master["Chr"].astype(str)
    vcf_all["Chr"] = vcf_all["Chr"].astype(str)

    merged = master.merge(
        vcf_all,
        on=["sample", "mode", "Chr", "Pos"],
        how="left",
        validate="many_to_one",
    )
    assert len(merged) == len(master), (
        f"Row count changed: master {len(master)} -> merged {len(merged)}"
    )

    # ---- Join-rate diagnostics ----
    join_stats = (
        merged.groupby(["sample", "mode"], observed=True)
        .agg(
            master_rows=("Pos", "size"),
            vcf_qual_hits=("vcf_QUAL", "count"),
            tp_label_match=("tp_label", "count"),
        )
        .assign(join_rate=lambda d: d["vcf_qual_hits"] / d["master_rows"])
        .reset_index()
    )
    join_path = LOG_DIR / "01_vcf_join_rate.tsv"
    join_stats.to_csv(join_path, sep="\t", index=False)
    print(f"[{time.strftime('%H:%M:%S')}] Wrote join stats -> {join_path}")
    print(join_stats.to_string(index=False))

    # ---- Sanity check: HCC1395 paired_full TP caller_AF vs ISM AF ----
    hcc = merged[
        (merged["sample"] == "HCC1395")
        & (merged["mode"] == "paired_full")
        & (merged["tp_label"] == 1)
    ].copy()
    hcc_valid = hcc.dropna(subset=["vcf_AF", "AF"])
    if len(hcc_valid) >= 100:
        r = hcc_valid["vcf_AF"].astype(float).corr(hcc_valid["AF"].astype(float))
        print(
            f"  HCC1395/paired_full TP rows with both AFs = {len(hcc_valid):,}; "
            f"Pearson(vcf_AF, AF) = {r:.4f}"
        )
    else:
        print("  HCC1395/paired_full TP too few rows for AF correlation check.")

    # ---- Write output ----
    out_path = OUT_DIR / "merged_with_vcf.tsv.gz"
    print(f"[{time.strftime('%H:%M:%S')}] Writing merged output -> {out_path}")
    merged.to_csv(out_path, sep="\t", index=False, compression="gzip")
    print(f"  final rows = {len(merged):,}; cols = {len(merged.columns)}")
    vcf_cols = [c for c in merged.columns if c.startswith("vcf_")]
    print(f"  vcf_* columns ({len(vcf_cols)}): {vcf_cols}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
